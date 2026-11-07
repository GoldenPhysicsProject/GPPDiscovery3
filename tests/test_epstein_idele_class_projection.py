"""Source-side exact class-character Fourier splitting over Q(sqrt(-5))."""
import unittest
from experiments.epstein_idele_class_projection import counts, report

class EpsteinClassTests(unittest.TestCase):
    def test_first_class_coefficients(self):
        a=counts(32,0);b=counts(32,1)
        self.assertEqual((a[2]//2,a[3]//2,a[6]//2),(0,0,2))
        self.assertEqual((b[2]//2,b[3]//2,b[6]//2),(1,2,0))
        self.assertEqual((a[1]//2,b[1]//2),(1,0))

    def test_fourier_channels_coprime_multiplicative_at_finite_limit(self):
        d=report(128)
        self.assertGreater(d["coprime_multiplicativity_pairs_tested"],300)
        self.assertEqual(d["nonmultiplicative_principal_witness"]["f(6)"],2)

    def test_even_representation_counts(self):
        for cls in (0,1):
            xs=counts(128,cls)
            self.assertTrue(all(v%2==0 for v in xs))

if __name__=="__main__":
    unittest.main()
