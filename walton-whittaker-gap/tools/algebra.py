"""Exact algebra behind the gap labels (Theorem B of FULL_PROOF.md); python integers / Fractions / Decimal only.

  * the characteristic polynomial of the size recurrence factors as (x - 1) f(x), f(x) = x^5 - 7x^4 + 15x^3 - 10x^2 - 3x + 2;
  * f is irreducible over Q: no root and no quadratic factor modulo 3 (exhaustive test of all monic polynomials of degree 1, 2);
  * Sturm's theorem: f has exactly three real roots; they are isolated by sign changes at rational points;
  * hence every root other than theta (the root in (3, 4)) has modulus < 3 < theta;
  * theta, 1/theta, 1 - 1/theta to 40 digits (exact bisection) and the polynomials they satisfy;
  * the recurrence reproduces the brute-force sizes of tests/sizes.py.
usage: python -m tools.algebra
"""
from decimal import Decimal, getcontext
from fractions import Fraction as Fr
from itertools import product

REC = [8, -22, 25, -7, -5, 2]                      # s_{k+6} = 8 s_{k+5} - 22 s_{k+4} + 25 s_{k+3} - 7 s_{k+2} - 5 s_{k+1} + 2 s_k
F = [2, -3, -10, 15, -7, 1]                        # f, coefficients of x^0 .. x^5


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def pev(p, x):
    r = 0
    for c in reversed(p):
        r = r * x + c
    return r


def trim(p):
    while p and p[-1] == 0:
        p = p[:-1]
    return p


def pmod(a, b, mod=None):
    a = [Fr(c) for c in a] if mod is None else [c % mod for c in a]
    b = trim([Fr(c) for c in b] if mod is None else [c % mod for c in b])
    a = trim(a)
    while len(a) >= len(b):
        if mod is None:
            f = a[-1] / b[-1]
        else:
            f = a[-1] * pow(b[-1], -1, mod) % mod
        sh = len(a) - len(b)
        for i, c in enumerate(b):
            a[sh + i] = a[sh + i] - f * c if mod is None else (a[sh + i] - f * c) % mod
        a = trim(a)
    return a


def sturm_count(p):
    """number of distinct real roots of p (Sturm's theorem, exact)"""
    seq = [[Fr(c) for c in p], [Fr(i * c) for i, c in enumerate(p)][1:]]
    while True:
        r = pmod(seq[-2], seq[-1])
        if not r:
            break
        seq.append([-c for c in r])

    def changes(signs):
        s = [x for x in signs if x != 0]
        return sum(1 for a, b in zip(s, s[1:]) if a != b)

    at_minus = [(1 if len(q) % 2 == 1 else -1) * (1 if q[-1] > 0 else -1) for q in seq]    # sign at -infinity
    at_plus = [1 if q[-1] > 0 else -1 for q in seq]                                       # sign at +infinity
    return changes(at_minus) - changes(at_plus)


def main():
    getcontext().prec = 60
    P6 = [-REC[5], -REC[4], -REC[3], -REC[2], -REC[1], -REC[0], 1]      # x^6 - 8x^5 + ... - 2, ascending coefficients
    assert P6 == pmul([-1, 1], F), "char poly = (x - 1) f(x)"
    print("x^6 - 8x^5 + 22x^4 - 25x^3 + 7x^2 + 5x - 2 = (x - 1)(x^5 - 7x^4 + 15x^3 - 10x^2 - 3x + 2): verified")

    # irreducibility mod 3
    f3 = [c % 3 for c in F]
    print("f mod 3 =", f3, "(ascending coefficients)")
    assert all(pev(f3, a) % 3 != 0 for a in range(3)), "root mod 3"
    for deg in (1, 2):
        for low in product(range(3), repeat=deg):
            g = list(low) + [1]
            assert trim(pmod(f3, g, 3)), f"f has the factor {g} modulo 3"
    print("f mod 3 has no monic factor of degree 1 or 2 (all 3 + 9 candidates tested): f is irreducible over F_3, hence over Q")

    # real roots
    n_real = sturm_count(F)
    assert n_real == 3
    print("Sturm: f has exactly", n_real, "real roots")
    pts = [Fr(-1, 2), Fr(-2, 5), Fr(39, 100), Fr(2, 5), Fr(3), Fr(4)]
    vals = [pev(F, x) for x in pts]
    print("f at -1/2, -2/5, 39/100, 2/5, 3, 4:", [str(v) for v in vals])
    assert vals[0] < 0 < vals[1] and vals[2] > 0 > vals[3] and vals[4] < 0 < vals[5]
    # product of the roots = -f(0) = -2: the two non-real roots z, conj(z) satisfy |z|^2 = 2 / (theta |r1| r2)
    bound = Fr(2) / (3 * Fr(2, 5) * Fr(39, 100))
    print(f"non-real roots: |z|^2 = 2/(theta |r1| r2) < 2/(3 * 2/5 * 39/100) = {bound} = {float(bound):.4f} < 9, so |z| < 3 < theta")
    assert bound < 9

    # theta by exact bisection
    lo, hi = Fr(3), Fr(4)
    for _ in range(140):
        mid = (lo + hi) / 2
        if pev(F, mid) < 0:
            lo = mid
        else:
            hi = mid
    th = (lo + hi) / 2

    def dec(q):
        return Decimal(q.numerator) / Decimal(q.denominator)

    print("theta        =", f"{dec(th):.40f}")
    print("1/theta      =", f"{dec(1 / th):.40f}")
    print("1 - 1/theta  =", f"{dec(1 - 1 / th):.40f}")

    # polynomials satisfied by y = 1/theta and l = 1 - 1/theta
    y_poly = list(reversed(F))                     # y^5 f(1/y)
    print("y = 1/theta is a root of", " + ".join(f"({c})y^{i}" for i, c in enumerate(y_poly)), "= 0")
    # l = 1 - y  <=>  y = 1 - l : substitute into y_poly (Horner with polynomial arithmetic in l)
    res = [0]
    for c in reversed(y_poly):
        res = pmul(res, [1, -1])
        res[0] += c
    res = trim(res)
    print("l = 1 - 1/theta is a root of", " + ".join(f"({c})l^{i}" for i, c in enumerate(res)), "= 0")
    assert abs(float(pev(res, 1 - 1 / th))) < 1e-12 and abs(float(pev(y_poly, 1 / th))) < 1e-12

    # the recurrence against the brute-force sizes of tests/sizes.py
    A = [1, 4, 16, 63, 245, 941, 3583, 13574, 51293, 193605]
    B = [2, 10, 42, 165, 632, 2396, 9049, 34135, 128737]
    for name, s in (('A', A), ('B', B)):
        pred = list(s[:6])
        for k in range(6, len(s)):
            pred.append(sum(c * pred[k - 1 - i] for i, c in enumerate(REC)))
        assert pred == s, name
        print(f"chain {name}: the recurrence started from s_0..s_5 reproduces the brute-force sizes s_6..s_{len(s) - 1}")
        print(f"   ratios s_(k-1)/s_k:", [f"{float(Fr(s[k - 1], s[k])):.5f}" for k in range(1, len(s))])


if __name__ == '__main__':
    main()
