"""Tiling rules, label windows and the 12-element symmetry group of the charged Walton-Whittaker tiling.

Lattice sites are pairs (a, b) = a*D0 + b*D1 with D_j = exp(i*pi*j/3).  A site carries a label (r, tau),
r in Z/6, tau = +-1.  Inflation of a site at p with label (r, tau) gives four children
    2p,            2p + D_{r+3}   labelled (r, tau)          (children 0, 1)
    2p + D_{r+2}                  labelled (r-1, tau)        (child 2)      -> children 0,1,2 form a triangle
    2p + D_{r+5}                  labelled (r+2, -tau)       (child 3, the "flipped" one).
"""
import math
import numpy as np

D = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]


def hexd(a, b):
    """lattice (hexagonal) distance of the offset (a, b)"""
    return max(abs(a), abs(b), abs(a + b))


def rot(p):
    """rotation by +60 degrees: D_j -> D_{j+1}"""
    return (-p[1], p[0] + p[1])


def rotk(p, k):
    for _ in range(k):
        p = rot(p)
    return p


def ball(R):
    return [(a, b) for a in range(-R, R + 1) for b in range(-R, R + 1) if hexd(a, b) <= R]


def euclid(p):
    a, b = p
    return math.sqrt((a + 0.5 * b) ** 2 + 0.75 * b * b)


def child_pos(o, r, j):
    """fine position of child j of the coarse site at offset o with label index r"""
    if j == 0:
        c = (0, 0)
    elif j == 1:
        c = D[(r + 3) % 6]
    elif j == 2:
        c = D[(r + 2) % 6]
    else:
        c = D[(r + 5) % 6]
    return (2 * o[0] + c[0], 2 * o[1] + c[1])


def child_label(r, t, j):
    return ((r % 6, t), (r % 6, t), ((r + 5) % 6, t), ((r + 2) % 6, -t))[j]


def build_patch(K, label=(0, 1)):
    """K-fold inflation of a single site at the origin: dict position -> label"""
    s = {(0, 0): label}
    for _ in range(K):
        new = {}
        for p, (r, t) in s.items():
            for j in range(4):
                q = child_pos(p, r, j)
                assert q not in new, "collision"
                new[q] = child_label(r, t, j)
        s = new
    return s


# ---------------------------------------------------------------------------------------------
# Windows: label content of the hexagonal ball of radius R around a site, as an int8 code vector.
# code 0 = no site; code = 1 + 4 r + 2 [tau = +1] + m  (m is a spare bit, always 0 here).
# ---------------------------------------------------------------------------------------------
def enc(r, t, m=0):
    return 1 + 4 * r + 2 * (1 if t == 1 else 0) + m


def dec(c):
    c -= 1
    return c >> 2, (1 if (c >> 1) & 1 else -1), c & 1


class Frame:
    """windows of radius R and the action of the group G = Z_6 (rotations) x Z_2 (species flip)"""

    def __init__(self, R):
        self.R = R
        self.offs = ball(R)
        self.idx = {o: i for i, o in enumerate(self.offs)}
        inv_rot = lambda p: (p[0] + p[1], -p[0])
        self.perm = np.array([self.idx[inv_rot(o)] for o in self.offs], dtype=np.int64)
        self.rotcode = np.zeros(25, dtype=np.int8)
        self.flipcode = np.zeros(25, dtype=np.int8)
        for c in range(1, 25):
            r, t, m = dec(c)
            self.rotcode[c] = enc((r + 1) % 6, t, m)
            self.flipcode[c] = enc(r, -t, m)

    def rotate(self, v):
        return self.rotcode[v[self.perm]]

    def flip(self, v):
        return self.flipcode[v]

    def canon(self, v):
        """(key, k, s): lexicographically least of the 12 images; k = number of 60-degree rotations used"""
        best = None
        cur = v
        for k in range(6):
            for s, a in ((0, cur), (1, self.flip(cur))):
                b = a.tobytes()
                if best is None or b < best[0]:
                    best = (b, k, s)
            cur = self.rotate(cur)
        return best

    def to_dict(self, v):
        return {o: dec(int(v[i])) for i, o in enumerate(self.offs) if v[i]}

    def from_dict(self, d):
        v = np.zeros(len(self.offs), dtype=np.int8)
        for o, (r, t, m) in d.items():
            if o in self.idx:
                v[self.idx[o]] = enc(r, t, m)
        return v

    def window_of(self, lab, p):
        """window of a labelled point set (dict position -> (r, t)) around position p"""
        v = np.zeros(len(self.offs), dtype=np.int8)
        for i, o in enumerate(self.offs):
            l = lab.get((p[0] + o[0], p[1] + o[1]))
            if l is not None:
                v[i] = enc(l[0], l[1])
        return v


def inflate(d):
    """d: dict coarse offset -> (r, t, m).  Returns dict fine position -> (r, t, j, parent offset).
    Asserts that no two children coincide (a collision would need two adjacent sites to claim the midpoint of their
    common edge); applied to every window of the closure this proves, by induction on the level, that the tiling is
    collision free: a collision at level K+1 would be seen inside the window of one of the two parents."""
    fine = {}
    for o, (r, t, m) in d.items():
        for j in range(4):
            lab = child_label(r, t, j)
            f = child_pos(o, r, j)
            assert f not in fine, "collision of two children"
            fine[f] = (lab[0], lab[1], j, o)
    return fine


def geometric_check(R):
    """Window of radius R of a child is determined by the radius-R window of the parent:
    every parent of a site within distance R of a child lies within distance R of the parent."""
    cj = [(0, 0)] + list(D)
    for c0 in cj:
        for o in ball(R + 3):
            if hexd(*o) <= R:
                continue
            for c in cj:
                if hexd(2 * o[0] + c[0] - c0[0], 2 * o[1] + c[1] - c0[1]) <= R:
                    return False
    return True
