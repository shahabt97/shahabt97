"""Finite set of inflation cells.

1. all-sites closure: every radius-R window of the tiling, obtained from the window of the origin by iterating
   the children map (a window of a child is a function of the window of its parent, Lemma 5 of FULL_PROOF.md);
2. seeds: windows centred at the member of a singleton cluster (chain A) or of a two-site cluster (chain B);
3. member closure: windows of the sites of C_k that are vertices of the triangle of a site z of C_{k-1} or flipped fourth
   children of neighbours of z, as a function of the window of z.  The least set containing the seeds and closed under
   this map is the set of cells (up to the symmetry group G).
"""
import math
import numpy as np
from .core import D, Frame, child_pos, euclid, inflate, build_patch

ZONE = 10.0   # fine sites farther than this (euclidean) from the inflated centre are not determined by a radius-7 window


# ---------------------------------------------------------------------------------------------
# 1. all windows of the tiling
# ---------------------------------------------------------------------------------------------
def origin_window(F, K=8, label=(0, 1)):
    patch = build_patch(K, label)
    w = {o: (patch[o][0], patch[o][1], 0) for o in F.offs if o in patch}
    patch2 = build_patch(K + 1, label)
    assert w == {o: (patch2[o][0], patch2[o][1], 0) for o in F.offs if o in patch2}, "origin window not stable"
    return w


def child_windows_plain(F, d):
    """windows of the four children of the centre of the window d"""
    fine = inflate(d)
    r0 = d[(0, 0)][0]
    out = []
    for j in range(4):
        cj = child_pos((0, 0), r0, j)
        w = {}
        for f, v in fine.items():
            o = (f[0] - cj[0], f[1] - cj[1])
            if o in F.idx:
                w[o] = (v[0], v[1], 0)
        out.append(w)
    return out


def all_sites_closure(F):
    d0 = origin_window(F)
    k0, _, _ = F.canon(F.from_dict(d0))
    seen = {k0: np.frombuffer(k0, dtype=np.int8).copy()}
    frontier = [seen[k0]]
    while frontier:
        nf = []
        for v in frontier:
            for w in child_windows_plain(F, F.to_dict(v)):
                k, _, _ = F.canon(F.from_dict(w))
                if k not in seen:
                    seen[k] = np.frombuffer(k, dtype=np.int8).copy()
                    nf.append(seen[k])
        frontier = nf
    return list(seen.values())


# ---------------------------------------------------------------------------------------------
# 2. seeds
# ---------------------------------------------------------------------------------------------
def same_species_nbrs(d, o):
    t = d[o][1]
    return [(o[0] + dl[0], o[1] + dl[1]) for dl in D
            if (o[0] + dl[0], o[1] + dl[1]) in d and d[(o[0] + dl[0], o[1] + dl[1])][1] == t]


def seed_windows(F, plain, interior_only=True):
    """('A', v): window centred at a singleton cluster; ('B', v): window centred at a site of a two-site cluster.
    interior_only: every lattice neighbour of every site of the seed cluster is a site of the tiling
    ("off the boundary"); the other seeds start chains next to holes of the tiling."""
    seeds = []
    z = (0, 0)
    for v in plain:
        d = F.to_dict(v)
        nb = same_species_nbrs(d, z)
        if len(nb) == 0:
            if interior_only and not all(dl in d for dl in D):
                continue
            seeds.append(('A', v))
        elif len(nb) == 1 and same_species_nbrs(d, nb[0]) == [z]:
            if interior_only and not all((s[0] + dl[0], s[1] + dl[1]) in d for s in (z, nb[0]) for dl in D):
                continue
            seeds.append(('B', v))
    return seeds


# ---------------------------------------------------------------------------------------------
# 3. member closure
# ---------------------------------------------------------------------------------------------
def component_from(fine, starts, zone=ZONE):
    """same-species connected component (nearest-neighbour graph) of the fine sites 'starts'"""
    comp = set(starts)
    fr = list(starts)
    while fr:
        nf = []
        for f in fr:
            sp = fine[f][1]
            for dl in D:
                g = (f[0] + dl[0], f[1] + dl[1])
                v = fine.get(g)
                if v is None or g in comp or v[1] != sp or euclid(g) > zone:
                    continue
                comp.add(g)
                nf.append(g)
        fr = nf
    return comp


def new_member_windows(F, d, bfs=False):
    """Windows of the possible new members associated with the centre z of d: the three triangle children of z and the
    flipped fourth children w_3 of the other-species neighbours w of z.  Every vertex of the cluster C_{k+1} that is
    not a vertex of a triangle T_z is such an extra w_3 with w adjacent to some z of C_k; the cell does not need to know
    whether w_3 really belongs to the cluster of z (that may be decided far outside the window), so all of them are
    returned (bfs=True returns only those found in the cluster of T_z inside the determined zone; used for counting)."""
    fine = inflate(d)
    z = (0, 0)
    r0, t0, _ = d[z]
    cands = [child_pos(z, r0, j) for j in range(3)]
    for q in D:
        if q in d and d[q][1] != t0:
            cands.append(child_pos(q, d[q][0], 3))
    if bfs:
        comp = component_from(fine, cands[:3])
        cands = [f for f in cands if f in comp]
    outs = []
    for f in cands:
        assert euclid(f) + F.R <= math.sqrt(3) * F.R - 1.0 + 1e-9, "output window not determined by input window"
        w = {}
        for o in F.offs:
            g = (f[0] + o[0], f[1] + o[1])
            v = fine.get(g)
            if v is not None:
                w[o] = (v[0], v[1], 0)
        outs.append(w)
    return outs


def member_closure(F, seeds, bfs=False):
    seen = {}
    frontier = []
    for _, v in seeds:
        k, _, _ = F.canon(v)
        if k not in seen:
            seen[k] = np.frombuffer(k, dtype=np.int8).copy()
            frontier.append(seen[k])
    while frontier:
        nf = []
        for v in frontier:
            for w in new_member_windows(F, F.to_dict(v), bfs):
                k, _, _ = F.canon(F.from_dict(w))
                if k not in seen:
                    seen[k] = np.frombuffer(k, dtype=np.int8).copy()
                    nf.append(seen[k])
        frontier = nf
    return list(seen.values())


def cells(R=7, interior_only=True):
    F = Frame(R)
    plain = all_sites_closure(F)
    seeds = seed_windows(F, plain, interior_only)
    return F, plain, seeds, member_closure(F, seeds)
