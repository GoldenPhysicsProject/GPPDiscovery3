"""Exact tests: no eigenvalue numerics, no gravity normalization assumed."""
import unittest
from math import log, pi, exp
from experiments.equator_fuzzy_sphere_car import equator_sl2, singlet_count, report

class EquatorFuzzySphereTests(unittest.TestCase):
    def test_exact_rotational_algebra(self):
        for N in range(1,13):
            r=equator_sl2(N)
            self.assertEqual(r["operator_dimension"],N*N)
            self.assertEqual(sum(2*l+1 for l in range(N)),N*N)

    def test_exact_SU2_singlet_counts(self):
        expected=[2,4,16,352,47040,32964544,108528014336,
                  1609374188955648]
        for N,val in enumerate(expected,1):
            self.assertEqual(singlet_count(N),val)

    def test_rigorous_lower_bound_finite_cases(self):
        c=2*exp(-1/24)/(3*pi**3)
        for N in range(1,13):
            d=singlet_count(N)
            self.assertGreaterEqual(d/2**(N*N),c/N**6)
            self.assertLessEqual(d,2**(N*N))

    def test_fermion_car_additional_not_matrix_entropy(self):
        d=report(3)
        self.assertEqual(d["rows"][2]["harmonics"],9)
        self.assertEqual(d["rows"][2]["Fock_dimension"],"512")
        self.assertEqual(d["rows"][2]["singlets"],"16")
        self.assertGreater(9*log(2),log(3))

if __name__=="__main__":
    unittest.main()
