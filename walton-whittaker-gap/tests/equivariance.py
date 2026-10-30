"""The children map on windows commutes with the group G = Z_6 x Z_2 (needed to work modulo G)."""
import random
import numpy as np
from ww.core import Frame, geometric_check
from ww.closure import all_sites_closure, child_windows_plain

F = Frame(7)
assert geometric_check(7)
plain = all_sites_closure(F)
random.seed(0)
bad = 0
for v in random.sample(plain, 200):
    ch = [F.from_dict(w) for w in child_windows_plain(F, F.to_dict(v))]
    vr = v
    for k in range(1, 6):
        vr = F.rotate(vr)
        for s in (0, 1):
            vv = F.flip(vr) if s else vr
            for a, w in zip(ch, child_windows_plain(F, F.to_dict(vv))):
                img = a
                for _ in range(k):
                    img = F.rotate(img)
                img = F.flip(img) if s else img
                bad += not np.array_equal(img, F.from_dict(w))
print("windows tested: 200 x 11 group elements x 4 children; equivariance violations:", bad)
assert bad == 0
