"""Exact projective-line/elliptic-curve zeta factor test."""
import unittest
from experiments.sphere_elliptic_zeta_bridge import (
    point_count_prime,point_count_quadratic,mirror_vs_height,verify
)
class SphereEllipticTests(unittest.TestCase):
    def test_elliptic_and_sphere_local_counts(self):
        for p in (3,5,7,11):
            n=point_count_prime(p); a=p+1-n
            self.assertEqual(point_count_quadratic(p),p*p+1-a*a+2*p)
    def test_mirror_does_not_imply_positive_height(self):
        r=mirror_vs_height()
        self.assertTrue(r["both_symplectic_dual"])
        self.assertFalse(r["bad_has_positive_invariant_height"])
    def test_full_exact_report(self):
        self.assertEqual(len(verify()["good_primes"]),4)
if __name__=="__main__":
    unittest.main()
