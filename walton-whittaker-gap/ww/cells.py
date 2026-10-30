"""An inflation cell: one inflation step applied to the radius-7 window of a site z of C_{k-1}.

It gives the labelled fine sites around 2z (inside the region determined by the window), their roles (child index, parent),
the triangles (children 0,1,2 of one parent), the extras (child 3) and the same-species nearest-neighbour graph.  The
two quadratic forms below are computed from it exactly.
"""
from fractions import Fraction as Fr
from .core import D, child_pos, inflate, euclid
from .closure import component_from, ZONE

# 20 * p5(x) with p5(x) = (2x^2+4x+1)^2 (x - 23/20) = 4x^5 + 57/5 x^4 + 8/5 x^3 - 15 x^2 - 41/5 x - 23/20
C20 = [-23, -164, -300, 32, 228, 80]          # coefficients of x^0 .. x^5


class Cell:
    def __init__(self, d):
        self.d = d
        self.fine = inflate(d)
        z = (0, 0)
        r0, t0, _ = d[z]
        self.tri = [child_pos(z, r0, j) for j in range(3)]
        self.comp = component_from(self.fine, self.tri)      # cluster of the triangle of z (inside the determined zone)
        # nearest-neighbour graph of ALL fine sites inside the zone: the graph ball around a site of a cluster is the
        # ball of this graph, whether or not the site is known (inside the zone) to belong to the cluster of T_z
        zone = [f for f in self.fine if euclid(f) <= ZONE]
        zs = set(zone)
        self.nbrs = {f: [(f[0] + dl[0], f[1] + dl[1]) for dl in D
                         if (f[0] + dl[0], f[1] + dl[1]) in zs and self.fine[(f[0] + dl[0], f[1] + dl[1])][1] == self.fine[f][1]]
                     for f in zone}
        self.tris = {}                      # parent offset -> fine sites of its triangle (cluster of T_z only)
        for f in self.comp:
            v = self.fine[f]
            if v[2] != 3:
                self.tris.setdefault(v[3], []).append(f)
        # fine sites whose rows of Q1 are checked in this cell: the triangle of z and EVERY flipped fourth child of an
        # other-species neighbour of z (whether it belongs to the cluster of z or to another one)
        self.owned = list(self.tri)
        for q in D:
            if q in d and d[q][1] != t0:
                self.owned.append(child_pos(q, d[q][0], 3))
        self._B = {}

    def dist_ball(self, x, rad):
        dist = {x: 0}
        fr = [x]
        for s in range(rad):
            nf = []
            for y in fr:
                assert euclid(y) <= ZONE - 1, "graph ball leaves the determined zone"
                for g in self.nbrs[y]:
                    if g not in dist:
                        dist[g] = s + 1
                        nf.append(g)
            fr = nf
        return dist

    def B(self, x):
        """B_x[v] = (A u_v)_x = number of neighbours of x inside the triangle of the coarse site v"""
        if x not in self._B:
            assert euclid(x) <= ZONE - 1, "site too close to the edge of the determined zone"
            out = {}
            for g in self.nbrs[x]:
                v = self.fine[g]
                if v[2] != 3:
                    out[v[3]] = out.get(v[3], 0) + 1
            self._B[x] = out
        return self._B[x]


# ------------------------------------------------------------------------------------------------
# lower edge:  Q1 = 9/8 I - A + 2 sum_v (A u_v)(A u_v)^T   (P = 2 sum_v u_v u_v^T)
# ------------------------------------------------------------------------------------------------
def q1_entry(cell, x, y, pw=2):
    val = Fr(0)
    if x == y:
        val += Fr(9, 8)
    if y in cell.nbrs[x]:
        val -= 1
    Bx, By = cell.B(x), cell.B(y)
    return val + pw * sum(c * By[v] for v, c in Bx.items() if v in By)


# ------------------------------------------------------------------------------------------------
# upper edge:  M = U^T p5(A) U  (U = triangle indicator vectors), rows computed as 20*M
# ------------------------------------------------------------------------------------------------
def coarse_data(cell):
    tri = {po: sites for po, sites in cell.tris.items() if len(sites) == 3}
    site2tri = {f: po for po, sites in tri.items() for f in sites}
    H = {po: set() for po in tri}
    for po, sites in tri.items():
        for f in sites:
            for g in cell.nbrs[f]:
                if g in site2tri and site2tri[g] != po:
                    H[po].add(site2tri[g])
    return tri, site2tri, H


def hball(H, x, rad, inside=None):
    """ball of radius rad in the triangle-adjacency graph H; inside(y) asserts that y is far enough from the zone edge"""
    dist = {x: 0}
    fr = [x]
    for s in range(rad):
        nf = []
        for y in fr:
            assert inside is None or inside(y), "H-ball leaves the determined zone"
            for g in H[y]:
                if g not in dist:
                    dist[g] = s + 1
                    nf.append(g)
        fr = nf
    return dist


def m_row(cell, tri, site2tri, po):
    """20 * M[po, w] for all coarse sites w: walks of length <= 5 starting in the triangle po"""
    a = [{f: 1 for f in tri[po]}]
    for _ in range(5):
        nxt = {}
        for f, v in a[-1].items():
            assert euclid(f) <= ZONE - 1, "walk leaves the determined zone"
            for g in cell.nbrs[f]:
                nxt[g] = nxt.get(g, 0) + v
        a.append(nxt)
    row = {}
    for i, ai in enumerate(a):
        for f, v in ai.items():
            w = site2tri.get(f)
            if w is not None:
                row[w] = row.get(w, 0) + C20[i] * v
    return row
