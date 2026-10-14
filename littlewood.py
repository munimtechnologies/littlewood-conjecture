"""Littlewood conjecture tools in pure Python (no dependencies).

Littlewood (about 1930) conjectured that for all real numbers x and y

    liminf_{q -> oo} q * ||q x|| * ||q y|| = 0,

where ||t|| is the distance from t to the nearest integer.

Real numbers are given exactly, by a function floor_scaled(N) = floor(N * x) computed with
integers only. So every continued fraction below is certified (not a floating-point guess),
and every Littlewood product comes with a proved error bound.

Commands (numbers: sqrt2, sqrt(7), cbrt2, 2^1/3, 4^1/3, golden, surd:a,b,d,c = (a + b sqrt d)/c):
    python littlewood.py cf X [N]         the first N partial quotients of X (default 30), certified
    python littlewood.py records X Y Q    every q <= Q at which q ||qX|| ||qY|| is a new minimum
    python littlewood.py value X Y q      a certified interval for q ||qX|| ||qY||
    python littlewood.py uniform X Y Q    T * min_{q <= T} ||qX|| ||qY|| for T = 1, 2, 4, ... <= Q
                                          (the uniform Littlewood quantity, disproved in 2026)
    python littlewood.py mixed X P Q      every q <= Q at which q |q|_P ||qX|| is a new minimum
                                          (the p-adic Littlewood conjecture of de Mathan-Teulie)

By MV2 of Munim, Inc. Made with AI agents working in MV2's name. MIT licence.
"""

import math
import re
import sys
from fractions import Fraction

DEFAULT_BITS = 256


# ---------------------------------------------------------------- exact real numbers

def iroot(n, k):
    """floor(n ** (1/k)) for integers n >= 0 and k >= 1, exactly (integer Newton iteration)."""
    if n < 0 or k < 1:
        raise ValueError("need n >= 0 and k >= 1")
    if k == 1 or n < 2:
        return n
    x = 1 << ((n.bit_length() + k - 1) // k)  # x > n ** (1/k)
    while True:
        y = ((k - 1) * x + n // x ** (k - 1)) // k
        if y >= x:
            return x
        x = y


def is_power(n, k):
    return n >= 0 and iroot(n, k) ** k == n


class Real:
    """A real number x given exactly by floor_scaled(N) = floor(N * x) for integers N >= 1."""

    def __init__(self, name, floor_scaled):
        self.name = name
        self.floor_scaled = floor_scaled
        self._approx = {}

    def approx(self, bits=DEFAULT_BITS):
        """A = floor(2^bits * x), so that A / 2^bits <= x < (A + 1) / 2^bits."""
        if bits not in self._approx:
            self._approx[bits] = self.floor_scaled(1 << bits)
        return self._approx[bits]

    def __float__(self):
        return self.approx(80) / (1 << 80)

    def __repr__(self):
        return self.name


def surd(a, b, d, c=1, name=None):
    """(a + b sqrt(d)) / c for integers with d > 0 not a square, b != 0 and c != 0."""
    if d <= 0 or is_power(d, 2):
        raise ValueError("d must be a positive non-square")
    if b == 0 or c == 0:
        raise ValueError("need b != 0 and c != 0")
    if c < 0:
        a, b, c = -a, -b, -c
    if name is None:
        name = "(%d + %d*sqrt(%d))/%d" % (a, b, d, c)

    def floor_scaled(N):
        # N x = (N a +- s) / c with s = sqrt(N^2 b^2 d) irrational, and floor(floor(t) / c) = floor(t / c)
        r = iroot(N * N * b * b * d, 2)
        top = N * a + r if b > 0 else N * a - r - 1
        return top // c

    return Real(name, floor_scaled)


def sqrt(d):
    return surd(0, 1, d, 1, name="sqrt(%d)" % d)


def root(p, k, q=1):
    """(p/q)^(1/k) for positive integers p, q and k >= 2, when it is irrational."""
    if p <= 0 or q <= 0 or k < 2:
        raise ValueError("need p, q > 0 and k >= 2")
    g = math.gcd(p, q)
    p, q = p // g, q // g
    if is_power(p, k) and is_power(q, k):
        raise ValueError("(%d/%d)^(1/%d) is rational" % (p, q, k))
    name = "%d^(1/%d)" % (p, k) if q == 1 else "(%d/%d)^(1/%d)" % (p, q, k)

    def floor_scaled(N):
        # floor(N x) = floor((N^k p / q)^(1/k)) = iroot(floor(N^k p / q), k)
        return iroot(N ** k * p // q, k)

    return Real(name, floor_scaled)


GOLDEN = surd(1, 1, 5, 2, name="golden ratio")


# ---------------------------------------------------------------- continued fractions

def cf_fraction(p, q):
    """Continued fraction of p/q (q > 0) by the floor algorithm."""
    if q <= 0:
        raise ValueError("need q > 0")
    out = []
    while q:
        a = p // q
        out.append(a)
        p, q = q, p - a * q
    return out


def continued_fraction(x, n, max_bits=1 << 16):
    """The first n partial quotients of x, certified.

    With A = floor(2^b x), x lies in [lo, hi] = [A / 2^b, (A + 1) / 2^b]. The reals whose floor
    algorithm starts a_0, ..., a_{j-1} form an interval, so every term shared by the expansions of
    lo and hi is a term of x. Precision doubles until n terms are shared."""
    bits = 64
    while bits <= max_bits:
        S = 1 << bits
        A = x.floor_scaled(S)
        common = []
        for s, t in zip(cf_fraction(A, S), cf_fraction(A + 1, S)):
            if s != t:
                break
            common.append(s)
        if len(common) >= n:
            return common[:n]
        bits *= 2
    raise ValueError("not enough precision for %d terms; is %s rational?" % (n, x))


def convergents(terms):
    """[(p_0, q_0), (p_1, q_1), ...] for the continued fraction [a_0; a_1, a_2, ...]."""
    out = []
    p_prev, q_prev, p, q = 0, 1, 1, 0
    for a in terms:
        p_prev, q_prev, p, q = p, q, a * p + p_prev, a * q + q_prev
        out.append((p, q))
    return out


# ---------------------------------------------------------------- distances and products

def dist(x, q, bits=DEFAULT_BITS):
    """||q x|| to within q / 2^bits."""
    S = 1 << bits
    r = q * x.approx(bits) % S
    return min(r, S - r) / S


def dist_bounds(x, q, bits=DEFAULT_BITS):
    """Fractions (lo, hi) with lo <= ||q x|| <= hi, for an integer q >= 1.

    q x lies in [q A / S, (q A + q) / S) with A = floor(S x), and ||.|| is 1-Lipschitz."""
    S = 1 << bits
    r = q * x.approx(bits) % S
    d = min(r, S - r)
    return Fraction(max(d - q, 0), S), Fraction(d + q, S)


def littlewood_value(x, y, q, bits=DEFAULT_BITS):
    """q ||q x|| ||q y|| as a float (relative error far below double precision for q << 2^100)."""
    S = 1 << bits
    rx = q * x.approx(bits) % S
    ry = q * y.approx(bits) % S
    return q * min(rx, S - rx) * min(ry, S - ry) / (S * S)


def littlewood_interval(x, y, q, bits=DEFAULT_BITS):
    """Fractions (lo, hi) with lo <= q ||q x|| ||q y|| <= hi."""
    xlo, xhi = dist_bounds(x, q, bits)
    ylo, yhi = dist_bounds(y, q, bits)
    return q * xlo * ylo, q * xhi * yhi


def _bits_for(Q):
    # each distance is then within q / 2^bits <= 2^-96 / Q of the truth
    return 2 * Q.bit_length() + 96


def littlewood_records(x, y, Q, bits=None):
    """Every q <= Q at which q ||q x|| ||q y|| is smaller than at every smaller q: [(q, value)]."""
    bits = bits or _bits_for(Q)
    S = 1 << bits
    ax, ay = x.approx(bits), y.approx(bits)
    rx = ry = 0
    best = None
    out = []
    for q in range(1, Q + 1):
        rx = (rx + ax) % S
        ry = (ry + ay) % S
        v = q * min(rx, S - rx) * min(ry, S - ry)
        if best is None or v < best:
            best = v
            out.append((q, v / (S * S)))
    return out


def uniform_profile(x, y, Q, bits=None):
    """[(T, T * min_{q <= T} ||q x|| ||q y||) for T = 1, 2, 4, ... <= Q].

    The uniform Littlewood conjecture asked that this tend to 0 for every pair; it is false
    (Schleischitz 2026). The classical conjecture only asks the same along some sequence of T."""
    bits = bits or _bits_for(Q)
    S = 1 << bits
    ax, ay = x.approx(bits), y.approx(bits)
    rx = ry = 0
    best = None
    out = []
    T = 1
    for q in range(1, Q + 1):
        rx = (rx + ax) % S
        ry = (ry + ay) % S
        v = min(rx, S - rx) * min(ry, S - ry)
        if best is None or v < best:
            best = v
        if q == T:
            out.append((T, T * best / (S * S)))
            T *= 2
    return out


def p_free_part(q, p):
    """q |q|_p: q with every factor p removed."""
    while q % p == 0:
        q //= p
    return q


def mixed_value(x, p, q, bits=DEFAULT_BITS):
    """q |q|_p ||q x|| as a float."""
    return p_free_part(q, p) * dist(x, q, bits)


def mixed_records(x, p, Q, bits=None):
    """Every q <= Q at which q |q|_p ||q x|| is a new minimum: [(q, value)]."""
    bits = bits or _bits_for(Q)
    S = 1 << bits
    ax = x.approx(bits)
    rx = 0
    best = None
    out = []
    for q in range(1, Q + 1):
        rx = (rx + ax) % S
        v = p_free_part(q, p) * min(rx, S - rx)
        if best is None or v < best:
            best = v
            out.append((q, v / S))
    return out


# ---------------------------------------------------------------- command line

def parse_real(s):
    t = s.strip().lower().replace(" ", "")
    if t in ("golden", "phi"):
        return GOLDEN
    m = re.fullmatch(r"sqrt\(?(\d+)\)?", t)
    if m:
        return sqrt(int(m.group(1)))
    m = re.fullmatch(r"cbrt\(?(\d+)\)?", t)
    if m:
        return root(int(m.group(1)), 3)
    m = re.fullmatch(r"(\d+)\^\(?(\d+)/(\d+)\)?", t)
    if m:
        n, a, k = (int(g) for g in m.groups())
        return root(n ** a, k)
    m = re.fullmatch(r"surd:(-?\d+),(-?\d+),(\d+),(-?\d+)", t)
    if m:
        return surd(*(int(g) for g in m.groups()))
    raise ValueError("cannot read %r: use sqrt2, cbrt2, 2^1/3, golden or surd:a,b,d,c" % s)


def main(argv):
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__)
        return 0
    cmd, args = argv[0], argv[1:]
    if cmd == "cf":
        x = parse_real(args[0])
        n = int(args[1]) if len(args) > 1 else 30
        terms = continued_fraction(x, n)
        print("%s = [%d; %s]" % (x, terms[0], ", ".join(str(a) for a in terms[1:])))
        if n > 1:
            print("largest of a_1 .. a_%d: %d" % (n - 1, max(terms[1:])))
        return 0
    if cmd == "records":
        x, y, Q = parse_real(args[0]), parse_real(args[1]), int(args[2])
        print("records of q ||q x|| ||q y||, x = %s, y = %s, q <= %d" % (x, y, Q))
        for q, v in littlewood_records(x, y, Q):
            print("%12d  %.12f" % (q, v))
        return 0
    if cmd == "value":
        x, y, q = parse_real(args[0]), parse_real(args[1]), int(args[2])
        lo, hi = littlewood_interval(x, y, q)
        print("%.17g <= %d ||%d x|| ||%d y|| <= %.17g   (x = %s, y = %s)" % (
            float(lo), q, q, q, float(hi), x, y))
        return 0
    if cmd == "uniform":
        x, y, Q = parse_real(args[0]), parse_real(args[1]), int(args[2])
        print("T * min_{q <= T} ||q x|| ||q y||, x = %s, y = %s" % (x, y))
        for T, v in uniform_profile(x, y, Q):
            print("%12d  %.12f" % (T, v))
        return 0
    if cmd == "mixed":
        x, p, Q = parse_real(args[0]), int(args[1]), int(args[2])
        print("records of q |q|_%d ||q x||, x = %s, q <= %d" % (p, x, Q))
        for q, v in mixed_records(x, p, Q):
            print("%12d  %.12f" % (q, v))
        return 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
