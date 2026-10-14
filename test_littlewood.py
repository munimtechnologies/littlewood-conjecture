"""Unit tests for littlewood.py. Run: python -m unittest test_littlewood (or python test_littlewood.py)."""

import contextlib
import io
import math
import unittest
from fractions import Fraction

import littlewood as lw

SQRT2 = lw.sqrt(2)
SQRT3 = lw.sqrt(3)
SQRT7 = lw.sqrt(7)
CBRT2 = lw.root(2, 3)

# OEIS A002945: continued fraction of 2^(1/3), terms a_0 .. a_40
CBRT2_CF = [1, 3, 1, 5, 1, 1, 4, 1, 1, 8, 1, 14, 1, 10, 2, 1, 4, 12, 2, 3, 2, 1, 3, 4, 1, 1, 2, 14,
            3, 12, 1, 15, 3, 1, 4, 534, 1, 1, 5, 1, 1]


class IntegerRoots(unittest.TestCase):
    def test_small(self):
        for k in range(1, 6):
            for n in range(0, 3000):
                r = lw.iroot(n, k)
                self.assertTrue(r ** k <= n < (r + 1) ** k, (n, k, r))

    def test_large(self):
        for n in (10 ** 60 + 12345, 2 ** 301 - 1, 3 ** 200, 3 ** 200 - 1):
            for k in (2, 3, 7):
                r = lw.iroot(n, k)
                self.assertTrue(r ** k <= n < (r + 1) ** k, (n, k))

    def test_is_power(self):
        self.assertTrue(lw.is_power(3 ** 200, 2))
        self.assertFalse(lw.is_power(3 ** 200 - 1, 2))
        self.assertTrue(lw.is_power(27, 3))
        self.assertFalse(lw.is_power(26, 3))


class ExactReals(unittest.TestCase):
    def test_floor_scaled_by_hand(self):
        self.assertEqual(SQRT2.floor_scaled(1000), 1414)
        self.assertEqual(lw.surd(1, -1, 2).floor_scaled(1000), -415)   # 1 - sqrt 2 = -0.414...
        self.assertEqual(lw.surd(1, 2, 2, 3).floor_scaled(3), 3)       # 1 + 2 sqrt 2 = 3.83...
        self.assertEqual(CBRT2.floor_scaled(10), 12)                   # 12.599...
        self.assertEqual(lw.GOLDEN.floor_scaled(1000), 1618)

    def test_floor_scaled_against_floats(self):
        cases = [
            (lw.surd(3, -2, 7, 5), (3 - 2 * math.sqrt(7)) / 5),
            (lw.surd(-4, 3, 11, -6), (-4 + 3 * math.sqrt(11)) / -6),
            (lw.root(5, 3, 2), 2.5 ** (1 / 3)),
            (lw.root(10, 5), 10 ** 0.2),
        ]
        for x, f in cases:
            for N in range(1, 300):
                self.assertEqual(x.floor_scaled(N), math.floor(N * f), (x, N))

    def test_float(self):
        self.assertAlmostEqual(float(SQRT2), math.sqrt(2), places=15)
        self.assertAlmostEqual(float(CBRT2), 2 ** (1 / 3), places=15)
        self.assertAlmostEqual(float(lw.GOLDEN), (1 + math.sqrt(5)) / 2, places=15)

    def test_rationals_rejected(self):
        for make in (lambda: lw.sqrt(4), lambda: lw.root(8, 3), lambda: lw.root(27, 3, 8),
                     lambda: lw.surd(1, 1, 9, 2), lambda: lw.surd(1, 0, 2), lambda: lw.sqrt(0)):
            with self.assertRaises(ValueError):
                make()


class ContinuedFractions(unittest.TestCase):
    def test_rational(self):
        self.assertEqual(lw.cf_fraction(415, 93), [4, 2, 6, 7])
        self.assertEqual(lw.convergents([4, 2, 6, 7])[-1], (415, 93))
        self.assertEqual(lw.cf_fraction(-7, 2), [-4, 2])
        self.assertEqual(lw.cf_fraction(5, 1), [5])

    def test_quadratic(self):
        self.assertEqual(lw.continued_fraction(SQRT2, 20), [1] + [2] * 19)
        self.assertEqual(lw.continued_fraction(SQRT3, 20), ([1] + [1, 2] * 10)[:20])
        self.assertEqual(lw.continued_fraction(SQRT7, 20), ([2] + [1, 1, 1, 4] * 5)[:20])
        self.assertEqual(lw.continued_fraction(lw.GOLDEN, 20), [1] * 20)
        self.assertEqual(lw.continued_fraction(lw.surd(1, -1, 2), 10), [-1, 1, 1] + [2] * 7)

    def test_cube_root_of_two(self):
        # matches OEIS A002945 through the partial quotient 534 (a_35) and beyond
        self.assertEqual(lw.continued_fraction(CBRT2, len(CBRT2_CF)), CBRT2_CF)

    def test_convergents_of_sqrt2(self):
        conv = lw.convergents(lw.continued_fraction(SQRT2, 11))
        self.assertEqual([q for p, q in conv], [1, 2, 5, 12, 29, 70, 169, 408, 985, 2378, 5741])
        self.assertEqual([p for p, q in conv], [1, 3, 7, 17, 41, 99, 239, 577, 1393, 3363, 8119])

    def test_legendre_bounds(self):
        # 1 / (q_{k+1} + q_k) < ||q_k x|| < 1 / q_{k+1} for k >= 1 (section 2 of the README)
        for x in (SQRT2, SQRT3, SQRT7, lw.GOLDEN, CBRT2, lw.surd(3, -2, 7, 5)):
            conv = lw.convergents(lw.continued_fraction(x, 16))
            for k in range(1, 15):
                q, q_next = conv[k][1], conv[k + 1][1]
                lo, hi = lw.dist_bounds(x, q)
                self.assertGreater(lo, Fraction(1, q_next + q), (x, k))
                self.assertLess(hi, Fraction(1, q_next), (x, k))

    def test_borel_hurwitz(self):
        # among any three consecutive convergents, one has q ||q x|| < 1/sqrt(5)
        for x in (SQRT2, SQRT3, SQRT7, lw.GOLDEN, CBRT2):
            conv = lw.convergents(lw.continued_fraction(x, 16))
            good = []
            for k in range(1, 16):
                q = conv[k][1]
                hi = q * lw.dist_bounds(x, q)[1]
                good.append(5 * hi * hi < 1)
            for k in range(len(good) - 2):
                self.assertTrue(any(good[k:k + 3]), (x, k))


class LittlewoodProducts(unittest.TestCase):
    def test_sqrt2_sqrt3_at_41(self):
        # 41 sqrt 2 = 57.98275..., 41 sqrt 3 = 71.01408...: 41 * 0.0172439 * 0.0140831 = 0.0099568
        v = lw.littlewood_value(SQRT2, SQRT3, 41)
        self.assertLess(abs(v - 0.0099568), 1e-6)
        lo, hi = lw.littlewood_interval(SQRT2, SQRT3, 41)
        self.assertLess(abs(float(lo) - v), 1e-16)
        self.assertLess(abs(float(hi) - v), 1e-16)
        self.assertLess(hi - lo, Fraction(1, 10 ** 30))
        self.assertLess(hi, Fraction(1, 26))
        records = lw.littlewood_records(SQRT2, SQRT3, 100)
        self.assertLessEqual(records[-1][1], v * (1 + 1e-12))

    def test_records_are_records(self):
        records = lw.littlewood_records(SQRT2, CBRT2, 3000)
        self.assertEqual(records[0][0], 1)
        for (q1, v1), (q2, v2) in zip(records, records[1:]):
            self.assertLess(q1, q2)
            self.assertLess(v2, v1)
        for q, v in records:
            self.assertLess(abs(v - lw.littlewood_value(SQRT2, CBRT2, q)), 1e-15)
        # no q between consecutive records beats the earlier record
        q_rec = [q for q, v in records]
        for q in range(1, 400):
            if q not in q_rec:
                best_before = min(v for r, v in records if r < q)
                self.assertGreaterEqual(lw.littlewood_value(SQRT2, CBRT2, q), best_before * (1 - 1e-12))

    def test_dependent_pair(self):
        # y = (1 + 2 sqrt 2) / 3, so 1, sqrt 2, y are linearly dependent; Lemma 3 with
        # q_k = 2378, q_{k+1} = 5741 and n = 3 q_k gives n ||n x|| ||n y|| <= 18 q_k / q_{k+1}^2 < 0.0013
        x, y = SQRT2, lw.surd(1, 2, 2, 3)
        self.assertLess(abs(lw.dist(y, 7134) - lw.dist(x, 4756)), 1e-12)
        lo, hi = lw.littlewood_interval(x, y, 7134)
        self.assertLess(hi, Fraction(13, 10000))
        records = lw.littlewood_records(x, y, 10000)
        self.assertLess(records[-1][1], 0.0013)

    def test_hurwitz_universal_bound(self):
        # m(x, y) <= 1 / (2 sqrt 5): some q has q ||q x|| < 1/sqrt 5, and ||q y|| <= 1/2
        for x, y in ((SQRT2, SQRT3), (lw.GOLDEN, CBRT2), (CBRT2, lw.root(4, 3)), (SQRT7, lw.GOLDEN)):
            best = lw.littlewood_records(x, y, 50)[-1][1]
            self.assertLess(best, 1 / (2 * math.sqrt(5)))

    def test_uniform_profile(self):
        prof = lw.uniform_profile(SQRT2, SQRT3, 1000)
        self.assertEqual([T for T, v in prof], [1, 2, 4, 8, 16, 32, 64, 128, 256, 512])
        first = lw.dist(SQRT2, 1) * lw.dist(SQRT3, 1)
        self.assertLess(abs(prof[0][1] - first), 1e-15)
        for T, v in prof:
            direct = T * min(lw.dist(SQRT2, q) * lw.dist(SQRT3, q) for q in range(1, T + 1))
            self.assertLess(abs(v - direct), 1e-12)

    def test_mixed(self):
        self.assertEqual(lw.p_free_part(12, 2), 3)
        self.assertEqual(lw.p_free_part(81, 3), 1)
        # 12 |12|_2 ||12 sqrt 2|| = 3 * 0.0294373 = 0.0883118
        self.assertLess(abs(lw.mixed_value(SQRT2, 2, 12) - 0.0883118), 1e-6)
        records = lw.mixed_records(SQRT2, 2, 2000)
        for (q1, v1), (q2, v2) in zip(records, records[1:]):
            self.assertLess(q1, q2)
            self.assertLess(v2, v1)
        for q, v in records:
            self.assertLess(abs(v - lw.mixed_value(SQRT2, 2, q)), 1e-15)


class CommandLine(unittest.TestCase):
    def test_parse(self):
        self.assertAlmostEqual(float(lw.parse_real("sqrt2")), math.sqrt(2), places=14)
        self.assertAlmostEqual(float(lw.parse_real("sqrt(3)")), math.sqrt(3), places=14)
        self.assertAlmostEqual(float(lw.parse_real("cbrt2")), 2 ** (1 / 3), places=14)
        self.assertAlmostEqual(float(lw.parse_real("2^1/3")), 2 ** (1 / 3), places=14)
        self.assertAlmostEqual(float(lw.parse_real("2^2/3")), 4 ** (1 / 3), places=14)
        self.assertAlmostEqual(float(lw.parse_real("golden")), (1 + math.sqrt(5)) / 2, places=14)
        self.assertAlmostEqual(float(lw.parse_real("surd:1,-1,2,1")), 1 - math.sqrt(2), places=14)
        with self.assertRaises(ValueError):
            lw.parse_real("pi")

    def test_commands_run(self):
        for argv in (["cf", "cbrt2", "12"], ["records", "sqrt2", "sqrt3", "100"],
                     ["value", "sqrt2", "sqrt3", "41"], ["uniform", "sqrt2", "sqrt3", "64"],
                     ["mixed", "sqrt2", "2", "100"]):
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                self.assertEqual(lw.main(argv), 0)
            self.assertTrue(out.getvalue().strip())
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            lw.main(["cf", "cbrt2", "12"])
        self.assertIn("2^(1/3) = [1; 3, 1, 5, 1, 1, 4, 1, 1, 8, 1, 14]", out.getvalue())


if __name__ == "__main__":
    unittest.main()
