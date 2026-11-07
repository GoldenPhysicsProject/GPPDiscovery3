"""Exact-source, floating-sign checks; not a directed interval certificate."""
import unittest
from experiments.prime_site_grade_signature_census import (
    site_atoms, site_matrix, leading_minors_3, inertia_from_ldl_3, run
)

class PrimeSiteSignatureTests(unittest.TestCase):
    def test_p3_power_support(self):
        self.assertEqual([a for a,_,_ in site_atoms(3,3)], [1,2,3,4,5])

    def test_ftheta_two_shift_source_normalization(self):
        z=site_matrix(3,3,-1,0)
        plus=site_matrix(3,3,-1,.25)
        minus=site_matrix(3,3,-1,-.25)
        both=site_matrix(3,3,-1,.25,include_conjugate=True)
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(both[i][j], plus[i][j]+minus[i][j],places=13)
                self.assertAlmostEqual(z[i][j], z[j][i],places=14)

    def test_p3_negative_parity_signature_flip(self):
        result=run(3,-1)
        z=result["matrices"]["zeta_weight_0"]
        f=result["matrices"]["Ftheta_combined"]
        self.assertEqual((z["positive"],z["negative"]),(2,1))
        self.assertEqual((f["positive"],f["negative"]),(3,0))
        self.assertLess(z["minors"][2],-0.00014)
        self.assertGreater(f["minors"][2],0.00010)
        self.assertEqual(result["num_powers"],5)

if __name__=="__main__":
    unittest.main()
