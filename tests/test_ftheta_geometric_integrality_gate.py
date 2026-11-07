"""Finite-field first-coefficient integrality gate."""
import unittest
from experiments.ftheta_geometric_integrality_gate import exact_gate
class FthetaGeometryTests(unittest.TestCase):
    def test_shifted_euler_count_not_integer(self):
        d=exact_gate()
        self.assertFalse(d["point_count_integral"])
        self.assertEqual(d["sphere_first_point_count"],3)
if __name__=="__main__":
    unittest.main()
