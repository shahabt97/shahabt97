"""Solve the piece SDP (numerically, SCS) and turn the solution into a rational certificate.

    maximise eps   s.t.   sum_c G_{tau(c)}[x,y] = F_xy - eps*[x=y]   (all constraints),   G_tau positive semidefinite.

The numerical solution is only a candidate: every piece is shifted by the smallest multiple of the identity that makes
it positive definite and rounded to a multiple of 2^-bits; ww.verify then checks everything exactly.

usage: python -m tools.solve_sdp constraints.pkl out.json.gz FORM [--tol 1e-4] [--scale 1] [--bits 30] [--margin 1e-6]
"""
import argparse
import gzip
import json
import pickle

import numpy as np
import scipy.sparse as sp
import cvxpy as cp


def solve(windows, sigs, tol, scale):
    sizes = [len(w) for w in windows]
    offs = [0]
    for n in sizes:
        offs.append(offs[-1] + n * n)
    merged = {}
    for lst, rhs, isdiag in sigs:
        nl = tuple(sorted((t, min(i, j), max(i, j)) for t, i, j in lst))
        merged[(nl, isdiag)] = rhs
    rows, cols, b, diag = [], [], [], []
    for r, ((nl, isdiag), rhs) in enumerate(merged.items()):
        for t, i, j in nl:
            rows.append(r)
            cols.append(offs[t] + j * sizes[t] + i)
        b.append(float(rhs) / scale)
        diag.append(1.0 if isdiag else 0.0)
    Aeq = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(merged), offs[-1]))
    G = [cp.Variable((n, n), PSD=True) for n in sizes]
    eps = cp.Variable()
    vec = cp.hstack([cp.vec(g, order='F') for g in G])
    prob = cp.Problem(cp.Maximize(eps), [Aeq @ vec + eps * np.array(diag) == np.array(b)])
    prob.solve(solver='SCS', eps=tol, max_iters=200000)
    return prob.status, float(eps.value) * scale, [g.value * scale for g in G]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('constraints')
    ap.add_argument('out')
    ap.add_argument('form', choices=['L', 'U'])
    ap.add_argument('--tol', type=float, default=1e-4)
    ap.add_argument('--scale', type=float, default=1.0)
    ap.add_argument('--bits', type=int, default=30)
    ap.add_argument('--margin', type=float, default=1e-6)
    a = ap.parse_args()
    data = pickle.load(open(a.constraints, 'rb'))
    status, eps, pieces = solve(data['windows'], data['sigs'], a.tol, a.scale)
    print(f"SDP status {status}, optimal slack eps* = {eps:.6g}")
    S = 1 << a.bits
    types = []
    for (j, key), win, G in zip(data['keys'], data['windows'], pieces):
        G = (G + G.T) / 2
        G = G + (max(0.0, -np.linalg.eigvalsh(G).min()) + a.margin) * np.eye(len(G))
        N = [[0] * len(G) for _ in G]
        for i in range(len(G)):
            for k in range(i, len(G)):
                N[i][k] = N[k][i] = int(round(float(G[i, k]) * S))
        types.append(dict(j=j, key=key.hex(), offs=[list(o) for o in win], N=N))
    with gzip.open(a.out, 'wt') as fh:
        json.dump(dict(form=a.form, Rt=data['Rt'], bits=a.bits, types=types), fh, separators=(',', ':'))
    print(f"{len(types)} pieces written to {a.out}")


if __name__ == '__main__':
    main()
